# Fargate Task Definitions Setup

This directory contains two Fargate task definitions for CosMx data processing:

1. **cosmx-process-slide** - Downloads slide from S3, runs stitch-images + read-targets sequentially, uploads results
2. **cosmx-metadata-generator** - Generates _metadata.csv with cell_id, cell_type, hex_color for all slides

## Well-Factored Architecture

### Common Patterns

Both task definitions share:
- **IAM Role**: `CosMxFargateProcessingRole` with S3 read/write permissions
- **Network Mode**: `awsvpc` (required for Fargate)
- **Logging**: CloudWatch Logs with consistent naming pattern
- **Region**: `us-west-2` (configurable)

### Resource Allocation Strategy

Tasks are sized based on computational requirements:

| Task | CPU | Memory | Ephemeral Storage | Rationale |
|------|-----|--------|-------------------|-----------|
| process-slide | 4 vCPU | 16 GB | 200 GB | Image + target processing (runs both sequentially) |
| metadata-generator | 1 vCPU | 2 GB | 20 GB (default) | Lightweight SQL processing |

### Storage Patterns

Both tasks use ephemeral storage and S3:
- **process-slide**: Downloads entire slide from S3 (CellStatsDir, AnalysisResults, RunSummary), runs stitch-images then read-targets sequentially, uploads results to `stitched/{slide_id}/`
- **metadata-generator**: Processes all slides, reads metadata files from S3 flatFiles, writes to `stitched/{slide_id}/_metadata.csv`

### Orchestration Pattern

For N slides, you'll launch **N + 1 tasks**:
- **N tasks** of `cosmx-process-slide` (one per slide, runs stitch-images + read-targets sequentially)
- **1 task** of `cosmx-metadata-generator` (processes all slides at once)

Each process-slide task:
1. Downloads all necessary slide data (CellStatsDir, RunSummary, AnalysisResults)
2. Auto-detects the AnalysisResults subdirectory (should be only one)
3. Runs stitch-images on the downloaded data
4. Runs read-targets on the same data (reuses downloads)
5. Uploads results to S3 stitched prefix

**Why sequential instead of parallel?**
- Can't filter S3 downloads by file extension (don't know what's needed upfront)
- Running sequentially allows reusing the same downloaded data for both operations
- Saves download time and ephemeral storage

## Prerequisites

Before registering these task definitions:

1. **IAM Roles**:
   ```bash
   # Task role (for S3 access)
   aws iam create-role \
     --role-name CosMxFargateProcessingRole \
     --assume-role-policy-document file://ecs-task-trust-policy.json

   aws iam put-role-policy \
     --role-name CosMxFargateProcessingRole \
     --policy-name CosmxS3Access \
     --policy-document file://cosmx-iam-policy.json

   # Execution role (for pulling images and logging)
   # Usually already exists, but if not:
   aws iam create-role \
     --role-name ecsTaskExecutionRole \
     --assume-role-policy-document file://ecs-task-trust-policy.json

   aws iam attach-role-policy \
     --role-name ecsTaskExecutionRole \
     --policy-arn arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy
   ```

2. **CloudWatch Log Groups**:
   ```bash
   aws logs create-log-group --log-group-name /ecs/cosmx-process-slide
   aws logs create-log-group --log-group-name /ecs/cosmx-metadata-generator
   ```

3. **Docker Images**:
   ```bash
   # Build and push headless image (for process-slide tasks)
   docker build --target headless -t ghcr.io/uw-brain-lab/cosmx-utilities:headless-latest .
   docker push ghcr.io/uw-brain-lab/cosmx-utilities:headless-latest

   # Build and push metadata-generator image (lightweight DuckDB-based)
   docker build --target metadata-generator -t ghcr.io/uw-brain-lab/cosmx-utilities:metadata-generator-latest .
   docker push ghcr.io/uw-brain-lab/cosmx-utilities:metadata-generator-latest
   ```

## Configuration

Infrastructure IDs (account, subnets, security group) are stored in `fargate/.env` (gitignored).

### First-time setup

```bash
cp fargate/.env.example fargate/.env
# Edit fargate/.env with your values
```

The `.env` file provides:

| Variable | Description |
|----------|-------------|
| `AWS_ACCOUNT_ID` | 12-digit AWS account ID |
| `AWS_REGION` | AWS region (default: `us-west-2`) |
| `ECS_CLUSTER` | ECS cluster name |
| `S3_BUCKET` | S3 bucket name for CosMx data |
| `FARGATE_SUBNETS` | Comma-separated subnet IDs |
| `FARGATE_SECURITY_GROUP` | Security group ID |

At runtime, the bucket and slide path are passed as CLI arguments to `process-slides.py`.

### Networking

Tasks run in the **cosmx VPC** on private subnets with a NAT Gateway for internet access.
`assignPublicIp=DISABLED` — the NAT Gateway handles outbound traffic (S3, GHCR, CloudWatch).

See `fargate/.env` for the specific subnet and security group IDs.

## Registration

Register both task definitions (substitutes `ACCOUNT_ID` from `.env`):

```bash
./fargate/register-task-defs.sh
```

## Running Tasks

### Option 1: AWS CLI

For a typical workflow with 3 slides, you'd launch 4 tasks total (3 process-slide + 1 metadata-generator):

```bash
source fargate/.env

BUCKET="${S3_BUCKET}"
BASE_PREFIX="CosMx-GBM/CosMx-GBM-segmentation-test-1.9.26/Bsbneurosegmentationtest1826_12_01_2026_12_26_12_282"
NETWORK_CONFIG="awsvpcConfiguration={subnets=[${FARGATE_SUBNETS}],securityGroups=[${FARGATE_SECURITY_GROUP}],assignPublicIp=DISABLED}"

# Launch process-slide for each slide (3 tasks)
# Each downloads entire slide, runs stitch-images + read-targets sequentially
for SLIDE_PATH in "${BASE_PREFIX}/slide-001" "${BASE_PREFIX}/slide-002" "${BASE_PREFIX}/slide-003"; do
  aws ecs run-task \
    --cluster "$ECS_CLUSTER" \
    --task-definition cosmx-process-slide \
    --launch-type FARGATE \
    --network-configuration "$NETWORK_CONFIG" \
    --overrides "{
      \"containerOverrides\": [{
        \"name\": \"process-slide\",
        \"command\": [\"bash\", \"/app/scripts/process-slide-wrapper.sh\", \"$BUCKET\", \"$SLIDE_PATH\"]
      }]
    }"
done

# Launch metadata generator once for all slides (1 task)
# Processes flatFiles for all slides and generates _metadata.csv with colors
aws ecs run-task \
  --cluster "$CLUSTER" \
  --task-definition cosmx-metadata-generator \
  --launch-type FARGATE \
  --network-configuration "$NETWORK_CONFIG" \
  --overrides "{
    \"containerOverrides\": [{
      \"name\": \"metadata-generator\",
      \"command\": [\"generate-metadata-csv.sh\", \"$BUCKET\", \"$BASE_PREFIX\"]
    }]
  }"
```

### Option 2: Step Functions (Recommended for Orchestration)

Create a Step Functions state machine to orchestrate tasks with proper dependencies.
Replace `<ECS_CLUSTER>`, `<FARGATE_SUBNETS>`, and `<FARGATE_SECURITY_GROUP>` with values from your `.env`:

```json
{
  "Comment": "CosMx Processing Pipeline - Process slides in parallel, then generate metadata",
  "StartAt": "ProcessSlides",
  "States": {
    "ProcessSlides": {
      "Type": "Map",
      "ItemsPath": "$.slides",
      "MaxConcurrency": 10,
      "Iterator": {
        "StartAt": "ProcessSingleSlide",
        "States": {
          "ProcessSingleSlide": {
            "Type": "Task",
            "Resource": "arn:aws:states:::ecs:runTask.sync",
            "Parameters": {
              "LaunchType": "FARGATE",
              "Cluster": "<ECS_CLUSTER>",
              "TaskDefinition": "cosmx-process-slide",
              "NetworkConfiguration": {
                "AwsvpcConfiguration": {
                  "Subnets": ["<FARGATE_SUBNETS>"],
                  "SecurityGroups": ["<FARGATE_SECURITY_GROUP>"],
                  "AssignPublicIp": "DISABLED"
                }
              },
              "Overrides": {
                "ContainerOverrides": [{
                  "Name": "process-slide",
                  "Command.$": "States.Array('bash', '/app/scripts/process-slide-wrapper.sh', $.bucket, $.slide_path)"
                }]
              }
            },
            "End": true
          }
        }
      },
      "Next": "GenerateMetadata"
    },
    "GenerateMetadata": {
      "Type": "Task",
      "Resource": "arn:aws:states:::ecs:runTask.sync",
      "Parameters": {
        "LaunchType": "FARGATE",
        "Cluster": "<ECS_CLUSTER>",
        "TaskDefinition": "cosmx-metadata-generator",
        "NetworkConfiguration": {
          "AwsvpcConfiguration": {
            "Subnets": ["<FARGATE_SUBNETS>"],
            "SecurityGroups": ["<FARGATE_SECURITY_GROUP>"],
            "AssignPublicIp": "DISABLED"
          }
        },
        "Overrides": {
          "ContainerOverrides": [{
            "Name": "metadata-generator",
            "Command.$": "States.Array('generate-metadata-csv.sh', $.bucket, $.base_prefix)"
          }]
        }
      },
      "End": true
    }
  }
}
```

**Input format for Step Functions:**
```json
{
  "bucket": "your-bucket-name",
  "base_prefix": "CosMx-GBM/experiment-123",
  "slides": [
    {
      "bucket": "your-bucket-name",
      "slide_path": "CosMx-GBM/experiment-123/slide-001"
    },
    {
      "bucket": "your-bucket-name",
      "slide_path": "CosMx-GBM/experiment-123/slide-002"
    }
  ]
}
```

## Maintenance

### Updating Task Definitions

When updating code:

1. Build and push new Docker image with version tag
2. Update task definition JSON with new image tag
3. Register new task definition revision
4. Update any automation to use new revision

### Cost Optimization

- **Spot Fargate**: For non-time-critical tasks, use `FARGATE_SPOT` capacity provider (up to 70% savings)
- **Right-sizing**: Monitor CloudWatch metrics and adjust CPU/memory if over/under-provisioned
- **Ephemeral storage**: Only process-slide needs 200GB; metadata-generator uses default 20GB

### Monitoring

Key metrics to watch in CloudWatch:
- Task CPU/Memory utilization (adjust if consistently under/over 80%)
- Task duration (optimize resource allocation)
- Failed task count (check logs for errors)
- S3 request metrics (download/upload patterns)
- Ephemeral storage usage (ensure 200GB is sufficient)

## Security Notes

1. **IAM Least Privilege**: The `CosMxFargateProcessingRole` only has S3 access to specific buckets/prefixes
2. **Network Isolation**: Tasks run on private subnets with NAT Gateway (no public IPs)
3. **No EFS**: All data in ephemeral storage (automatically encrypted, deleted after task)
4. **Secrets**: Use AWS Secrets Manager for any sensitive configuration (not needed currently)

## Troubleshooting

### Task fails to start
- Check execution role has permission to pull Docker image from GHCR
- Verify log group exists in CloudWatch Logs
- Ensure security group allows outbound HTTPS (port 443 for ECR and S3)

### Task fails during execution
- Check CloudWatch logs for the specific task (use stream prefix to find logs)
- Verify IAM task role has required S3 permissions (GetObject, PutObject, ListBucket)
- Check S3 bucket/prefix exists and is accessible

### Out of memory errors
- Increase memory allocation in task definition (currently 16GB for process-slide)
- Monitor actual memory usage in CloudWatch to right-size

### Out of ephemeral storage
- Increase ephemeralStorage.sizeInGiB (currently 200GB, max is 200GB)
- Check if intermediate files are being cleaned up properly
- Consider processing fewer FOVs per slide if data is too large

### AnalysisResults subdirectory not found
- Verify the subdirectory exists in S3 under AnalysisResults/
- Check CloudWatch logs for the detected subdirectory name
- Ensure aws cli has permissions to list S3 objects

## References

- [AWS Fargate Task Definitions](https://docs.aws.amazon.com/AmazonECS/latest/developerguide/task_definitions.html)
- [ECS Task IAM Roles](https://docs.aws.amazon.com/AmazonECS/latest/developerguide/task-iam-roles.html)
- [Fargate Ephemeral Storage](https://docs.aws.amazon.com/AmazonECS/latest/developerguide/fargate-task-storage.html)
